% Calculate inertial parameters regressor of gravitation load for
% S5RPRPR9
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
% Datum: 2019-12-31 18:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRPR9_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR9_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR9_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR9_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t25 = sin(qJ(3));
t23 = qJ(1) + pkin(8);
t17 = sin(t23);
t18 = cos(t23);
t9 = g(1) * t18 + g(2) * t17;
t50 = t9 * t25;
t19 = t25 * qJ(4);
t28 = cos(qJ(3));
t38 = t28 * pkin(3) + t19;
t2 = g(3) * t25 + t9 * t28;
t48 = pkin(3) * t25;
t47 = g(1) * t17;
t43 = g(3) * t28;
t42 = t28 * pkin(7);
t41 = t18 * t28;
t24 = sin(qJ(5));
t40 = t24 * t25;
t27 = cos(qJ(5));
t39 = t25 * t27;
t37 = qJ(4) * t28;
t29 = cos(qJ(1));
t36 = t29 * pkin(1) + t18 * pkin(2) + t17 * pkin(6);
t26 = sin(qJ(1));
t35 = -t26 * pkin(1) + t18 * pkin(6);
t34 = pkin(3) * t41 + t18 * t19 + t36;
t33 = -g(2) * t18 + t47;
t32 = g(1) * t26 - g(2) * t29;
t31 = -pkin(2) - t38;
t12 = t18 * t37;
t10 = t17 * t37;
t8 = t33 * t28;
t7 = t33 * t25;
t6 = -t17 * t40 + t18 * t27;
t5 = t17 * t39 + t18 * t24;
t4 = t17 * t27 + t18 * t40;
t3 = -t17 * t24 + t18 * t39;
t1 = -t43 + t50;
t11 = [0, 0, 0, 0, 0, 0, t32, g(1) * t29 + g(2) * t26, 0, 0, 0, 0, 0, 0, 0, 0, t33, t9, 0, t32 * pkin(1), 0, 0, 0, 0, 0, 0, t8, -t7, -t9, -g(1) * (-t17 * pkin(2) + t35) - g(2) * t36, 0, 0, 0, 0, 0, 0, -t9, -t8, t7, -g(1) * t35 - g(2) * t34 - t31 * t47, 0, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t4, g(1) * t5 - g(2) * t3, t8, -g(1) * (t18 * pkin(4) + t35) - g(2) * (pkin(7) * t41 + t34) + (-g(1) * (t31 - t42) - g(2) * pkin(4)) * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, -t2, -g(1) * (-t18 * t48 + t12) - g(2) * (-t17 * t48 + t10) - g(3) * t38, 0, 0, 0, 0, 0, 0, -t2 * t24, -t2 * t27, t1, -g(1) * t12 - g(2) * t10 - g(3) * (t38 + t42) + (pkin(3) + pkin(7)) * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t3 - g(2) * t5 + t27 * t43, g(1) * t4 - g(2) * t6 - t24 * t43, 0, 0;];
taug_reg = t11;
