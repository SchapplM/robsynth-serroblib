% Calculate inertial parameters regressor of gravitation load for
% S5RRPPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPPR11_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR11_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR11_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR11_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t27 = sin(qJ(2));
t28 = sin(qJ(1));
t30 = cos(qJ(1));
t9 = g(1) * t30 + g(2) * t28;
t60 = t9 * t27;
t17 = t27 * qJ(3);
t29 = cos(qJ(2));
t40 = t29 * pkin(2) + t17;
t6 = g(3) * t27 + t29 * t9;
t59 = pkin(2) * t27;
t24 = sin(pkin(8));
t58 = pkin(4) * t24;
t57 = g(1) * t28;
t53 = g(3) * t29;
t23 = pkin(8) + qJ(5);
t15 = sin(t23);
t51 = t28 * t15;
t16 = cos(t23);
t50 = t28 * t16;
t49 = t28 * t24;
t25 = cos(pkin(8));
t48 = t28 * t25;
t26 = -pkin(7) - qJ(4);
t47 = t29 * t26;
t46 = t29 * t30;
t45 = t30 * t15;
t44 = t30 * t16;
t43 = t30 * t24;
t42 = t30 * t25;
t39 = t30 * pkin(1) + t28 * pkin(6);
t38 = qJ(3) * t29;
t37 = t29 * qJ(4);
t36 = t29 * t58;
t35 = t27 * t43;
t34 = pkin(2) * t46 + t30 * t17 + t39;
t33 = -g(2) * t30 + t57;
t32 = -pkin(1) - t40;
t20 = t30 * pkin(6);
t14 = pkin(4) * t25 + pkin(3);
t12 = t30 * t38;
t10 = t28 * t38;
t8 = t33 * t29;
t7 = t33 * t27;
t5 = -t53 + t60;
t4 = -t27 * t51 + t44;
t3 = t27 * t50 + t45;
t2 = t27 * t45 + t50;
t1 = t27 * t44 - t51;
t11 = [0, 0, 0, 0, 0, 0, t33, t9, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t7, -t9, -g(1) * (-t28 * pkin(1) + t20) - g(2) * t39, 0, 0, 0, 0, 0, 0, -t9, -t8, t7, -g(1) * t20 - g(2) * t34 - t32 * t57, 0, 0, 0, 0, 0, 0, -g(1) * (-t27 * t49 + t42) - g(2) * (t35 + t48), -g(1) * (-t27 * t48 - t43) - g(2) * (t27 * t42 - t49), t8, -g(1) * (t30 * pkin(3) + t20) - g(2) * (t30 * t37 + t34) + (-g(1) * (t32 - t37) - g(2) * pkin(3)) * t28, 0, 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t2, g(1) * t3 - g(2) * t1, t8, -g(1) * (t30 * t14 + t20) - g(2) * (pkin(4) * t35 - t26 * t46 + t34) + (-g(1) * (-t27 * t58 + t32 + t47) - g(2) * t14) * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, -t6, -g(1) * (-t30 * t59 + t12) - g(2) * (-t28 * t59 + t10) - g(3) * t40, 0, 0, 0, 0, 0, 0, -t6 * t24, -t6 * t25, t5, -g(1) * t12 - g(2) * t10 - g(3) * (t37 + t40) + (pkin(2) + qJ(4)) * t60, 0, 0, 0, 0, 0, 0, -t6 * t15, -t6 * t16, t5, -g(1) * (t30 * t36 + t12) - g(2) * (t28 * t36 + t10) - g(3) * (t40 - t47) + (-g(3) * t58 + t9 * (pkin(2) - t26)) * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 - g(2) * t3 + t16 * t53, g(1) * t2 - g(2) * t4 - t15 * t53, 0, 0;];
taug_reg = t11;
