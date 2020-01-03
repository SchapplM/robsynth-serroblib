% Calculate inertial parameters regressor of gravitation load for
% S4PRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% taug_reg [4x(4*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4PRRR7_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR7_gravloadJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR7_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4PRRR7_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t20 = sin(pkin(4));
t46 = g(3) * t20;
t23 = sin(qJ(2));
t45 = t20 * t23;
t25 = cos(qJ(3));
t44 = t20 * t25;
t26 = cos(qJ(2));
t43 = t20 * t26;
t21 = sin(qJ(4));
t42 = t21 * t25;
t24 = cos(qJ(4));
t41 = t24 * t25;
t40 = t25 * t26;
t39 = pkin(2) * t43 + pkin(6) * t45;
t38 = cos(pkin(4));
t37 = cos(pkin(8));
t19 = sin(pkin(8));
t31 = t38 * t37;
t8 = t19 * t23 - t26 * t31;
t9 = t19 * t26 + t23 * t31;
t36 = -t8 * pkin(2) + t9 * pkin(6);
t34 = t19 * t38;
t10 = t37 * t23 + t26 * t34;
t11 = -t23 * t34 + t37 * t26;
t35 = -t10 * pkin(2) + t11 * pkin(6);
t33 = t20 * t37;
t22 = sin(qJ(3));
t32 = pkin(3) * t25 + pkin(7) * t22;
t12 = -t22 * t45 + t38 * t25;
t2 = -t9 * t22 - t25 * t33;
t4 = -t11 * t22 + t19 * t44;
t30 = g(1) * t4 + g(2) * t2 + g(3) * t12;
t13 = t38 * t22 + t23 * t44;
t3 = -t22 * t33 + t9 * t25;
t5 = t19 * t20 * t22 + t11 * t25;
t29 = g(1) * t5 + g(2) * t3 + g(3) * t13;
t28 = -g(1) * t10 - g(2) * t8 + g(3) * t43;
t27 = g(1) * t11 + g(2) * t9 + g(3) * t45;
t1 = t28 * t22;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, t27, 0, 0, 0, 0, 0, 0, 0, 0, -t28 * t25, t1, -t27, -g(1) * t35 - g(2) * t36 - g(3) * t39, 0, 0, 0, 0, 0, 0, -g(1) * (-t10 * t41 + t11 * t21) - g(2) * (t9 * t21 - t8 * t41) - (t21 * t23 + t24 * t40) * t46, -g(1) * (t10 * t42 + t11 * t24) - g(2) * (t9 * t24 + t8 * t42) - (-t21 * t40 + t23 * t24) * t46, -t1, -g(1) * (-t32 * t10 + t35) - g(2) * (-t32 * t8 + t36) - g(3) * (t32 * t43 + t39); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, t29, 0, 0, 0, 0, 0, 0, 0, 0, -t30 * t24, t30 * t21, -t29, -g(1) * (t4 * pkin(3) + t5 * pkin(7)) - g(2) * (t2 * pkin(3) + t3 * pkin(7)) - g(3) * (t12 * pkin(3) + t13 * pkin(7)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t10 * t24 - t5 * t21) - g(2) * (-t3 * t21 + t8 * t24) - g(3) * (-t13 * t21 - t24 * t43), -g(1) * (-t10 * t21 - t5 * t24) - g(2) * (-t8 * t21 - t3 * t24) - g(3) * (-t13 * t24 + t21 * t43), 0, 0;];
taug_reg = t6;
