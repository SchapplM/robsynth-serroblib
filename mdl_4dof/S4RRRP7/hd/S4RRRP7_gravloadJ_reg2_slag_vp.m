% Calculate inertial parameters regressor of gravitation load for
% S4RRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% taug_reg [4x(4*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RRRP7_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP7_gravloadJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP7_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP7_gravloadJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t25 = sin(qJ(1));
t28 = cos(qJ(1));
t11 = g(1) * t28 + g(2) * t25;
t24 = sin(qJ(2));
t48 = t11 * t24;
t27 = cos(qJ(2));
t38 = t27 * pkin(2) + t24 * pkin(6);
t47 = g(1) * t25;
t44 = g(3) * t24;
t43 = t24 * t28;
t42 = t25 * t27;
t41 = t27 * t28;
t23 = sin(qJ(3));
t40 = t28 * t23;
t26 = cos(qJ(3));
t39 = t28 * t26;
t37 = t28 * pkin(1) + t25 * pkin(5);
t36 = pkin(2) * t41 + pkin(6) * t43 + t37;
t6 = t23 * t42 + t39;
t8 = -t25 * t26 + t27 * t40;
t35 = g(1) * t6 - g(2) * t8;
t34 = -g(2) * t28 + t47;
t33 = pkin(3) * t26 + qJ(4) * t23;
t1 = g(1) * t8 + g(2) * t6 + t23 * t44;
t7 = t26 * t42 - t40;
t9 = t25 * t23 + t27 * t39;
t31 = g(1) * t9 + g(2) * t7 + t26 * t44;
t30 = (-pkin(1) - t38) * t47;
t29 = -g(3) * t27 + t48;
t21 = t28 * pkin(5);
t14 = pkin(6) * t41;
t12 = pkin(6) * t42;
t10 = t34 * t24;
t5 = t11 * t27 + t44;
t4 = t29 * t26;
t3 = t29 * t23;
t2 = g(1) * t7 - g(2) * t9;
t13 = [0, 0, 0, 0, 0, 0, t34, t11, 0, 0, 0, 0, 0, 0, 0, 0, t34 * t27, -t10, -t11, -g(1) * (-t25 * pkin(1) + t21) - g(2) * t37, 0, 0, 0, 0, 0, 0, t2, -t35, t10, -g(1) * t21 - g(2) * t36 - t30, 0, 0, 0, 0, 0, 0, t2, t10, t35, -g(1) * (-t7 * pkin(3) - t6 * qJ(4) + t21) - g(2) * (t9 * pkin(3) + t8 * qJ(4) + t36) - t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, t5, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t3, -t5, -g(1) * (-pkin(2) * t43 + t14) - g(2) * (-t25 * t24 * pkin(2) + t12) - g(3) * t38, 0, 0, 0, 0, 0, 0, t4, -t5, t3, -g(1) * t14 - g(2) * t12 - g(3) * (t33 * t27 + t38) + (pkin(2) + t33) * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t31, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t31, -g(1) * (-t8 * pkin(3) + t9 * qJ(4)) - g(2) * (-t6 * pkin(3) + t7 * qJ(4)) - (-pkin(3) * t23 + qJ(4) * t26) * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t13;
