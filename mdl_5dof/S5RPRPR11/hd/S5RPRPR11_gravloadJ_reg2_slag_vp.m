% Calculate inertial parameters regressor of gravitation load for
% S5RPRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRPR11_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR11_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR11_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR11_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t31 = sin(qJ(1));
t33 = cos(qJ(1));
t17 = g(1) * t33 + g(2) * t31;
t26 = pkin(8) + qJ(3);
t23 = sin(t26);
t49 = t17 * t23;
t19 = t23 * qJ(4);
t24 = cos(t26);
t40 = t24 * pkin(3) + t19;
t45 = t24 * pkin(4);
t29 = -pkin(6) - qJ(2);
t44 = pkin(7) + t29;
t43 = t23 * t33;
t42 = t24 * t33;
t41 = t33 * t29;
t39 = qJ(4) * t24;
t28 = cos(pkin(8));
t22 = t28 * pkin(2) + pkin(1);
t15 = t33 * t22;
t38 = g(2) * (pkin(3) * t42 + t33 * t19 + t15);
t16 = g(1) * t31 - g(2) * t33;
t30 = sin(qJ(5));
t32 = cos(qJ(5));
t10 = t23 * t32 - t24 * t30;
t37 = t23 * t30 + t24 * t32;
t36 = -t22 - t40;
t3 = t10 * t31;
t5 = t30 * t42 - t32 * t43;
t35 = g(1) * t5 - g(2) * t3 + g(3) * t37;
t4 = t37 * t31;
t6 = t37 * t33;
t34 = g(1) * t6 + g(2) * t4 + g(3) * t10;
t14 = t33 * t39;
t12 = t31 * t39;
t8 = t16 * t24;
t7 = t16 * t23;
t2 = g(3) * t23 + t17 * t24;
t1 = -g(3) * t24 + t49;
t9 = [0, 0, 0, 0, 0, 0, t16, t17, 0, 0, 0, 0, 0, 0, 0, 0, t16 * t28, -t16 * sin(pkin(8)), -t17, -g(1) * (-t31 * pkin(1) + t33 * qJ(2)) - g(2) * (t33 * pkin(1) + t31 * qJ(2)), 0, 0, 0, 0, 0, 0, t8, -t7, -t17, -g(1) * (-t31 * t22 - t41) - g(2) * (-t31 * t29 + t15), 0, 0, 0, 0, 0, 0, t8, -t17, t7, g(1) * t41 - t38 + (-g(1) * t36 + g(2) * t29) * t31, 0, 0, 0, 0, 0, 0, g(1) * t4 - g(2) * t6, g(1) * t3 + g(2) * t5, t17, -t38 + (g(1) * t44 - g(2) * t45) * t33 + (-g(1) * (t36 - t45) + g(2) * t44) * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t2, -g(1) * (-pkin(3) * t43 + t14) - g(2) * (-t31 * t23 * pkin(3) + t12) - g(3) * t40, 0, 0, 0, 0, 0, 0, -t35, -t34, 0, -g(1) * t14 - g(2) * t12 - g(3) * (t40 + t45) + (pkin(3) + pkin(4)) * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, t34, 0, 0;];
taug_reg = t9;
