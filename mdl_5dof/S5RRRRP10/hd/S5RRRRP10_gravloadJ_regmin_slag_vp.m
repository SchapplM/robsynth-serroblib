% Calculate minimal parameter regressor of gravitation load for
% S5RRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
% 
% Output:
% taug_reg [5x26]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRRP10_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP10_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP10_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP10_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t26 = sin(qJ(2));
t27 = sin(qJ(1));
t30 = cos(qJ(2));
t41 = cos(pkin(5));
t51 = cos(qJ(1));
t36 = t41 * t51;
t12 = t27 * t26 - t30 * t36;
t24 = sin(qJ(4));
t28 = cos(qJ(4));
t13 = t26 * t36 + t27 * t30;
t25 = sin(qJ(3));
t29 = cos(qJ(3));
t22 = sin(pkin(5));
t39 = t22 * t51;
t5 = t13 * t29 - t25 * t39;
t58 = -t12 * t28 + t5 * t24;
t57 = t12 * t24 + t5 * t28;
t38 = t27 * t41;
t14 = t51 * t26 + t30 * t38;
t15 = -t26 * t38 + t51 * t30;
t47 = t22 * t27;
t9 = t15 * t29 + t25 * t47;
t1 = t14 * t28 - t9 * t24;
t46 = t22 * t29;
t11 = t41 * t25 + t26 * t46;
t45 = t22 * t30;
t56 = g(2) * t58 - g(3) * (-t11 * t24 - t28 * t45) - g(1) * t1;
t54 = g(2) * t12;
t53 = g(2) * t13;
t52 = g(3) * t22;
t48 = t22 * t26;
t44 = t24 * t29;
t43 = t28 * t29;
t42 = t29 * t30;
t40 = pkin(4) * t24 + pkin(8);
t4 = t13 * t25 + t29 * t39;
t8 = t15 * t25 - t27 * t46;
t37 = -g(1) * t4 + g(2) * t8;
t21 = t28 * pkin(4) + pkin(3);
t23 = -qJ(5) - pkin(9);
t35 = t21 * t29 - t23 * t25 + pkin(2);
t10 = t25 * t48 - t41 * t29;
t34 = g(1) * t8 + g(2) * t4 + g(3) * t10;
t33 = g(1) * t9 + g(2) * t5 + g(3) * t11;
t32 = -g(1) * t14 + g(3) * t45 - t54;
t3 = t32 * t25;
t2 = t14 * t24 + t9 * t28;
t6 = [0, g(1) * t27 - g(2) * t51, g(1) * t51 + g(2) * t27, 0, 0, 0, 0, 0, g(1) * t13 - g(2) * t15, -g(1) * t12 + g(2) * t14, 0, 0, 0, 0, 0, g(1) * t5 - g(2) * t9, t37, 0, 0, 0, 0, 0, g(1) * t57 - g(2) * t2, -g(1) * t58 - g(2) * t1, -t37, -g(1) * (-t27 * pkin(1) - t13 * pkin(2) + pkin(7) * t39 - t40 * t12 - t21 * t5 + t23 * t4) - g(2) * (t51 * pkin(1) + t15 * pkin(2) + pkin(7) * t47 + t40 * t14 + t9 * t21 - t8 * t23); 0, 0, 0, 0, 0, 0, 0, 0, -t32, g(1) * t15 + g(3) * t48 + t53, 0, 0, 0, 0, 0, -t32 * t29, t3, 0, 0, 0, 0, 0, -g(1) * (-t14 * t43 + t15 * t24) - g(2) * (-t12 * t43 + t13 * t24) - (t24 * t26 + t28 * t42) * t52, -g(1) * (t14 * t44 + t15 * t28) - g(2) * (t12 * t44 + t13 * t28) - (-t24 * t42 + t26 * t28) * t52, -t3, -g(1) * (-t35 * t14 + t40 * t15) - t40 * t53 + t35 * t54 - (t40 * t26 + t35 * t30) * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, t33, 0, 0, 0, 0, 0, t34 * t28, -t34 * t24, -t33, -g(1) * (-t8 * t21 - t9 * t23) - g(2) * (-t4 * t21 - t5 * t23) - g(3) * (-t10 * t21 - t11 * t23); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, g(1) * t2 + g(2) * t57 - g(3) * (-t11 * t28 + t24 * t45), 0, t56 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34;];
taug_reg = t6;
