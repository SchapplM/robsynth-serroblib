% Calculate inertial parameters regressor of gravitation load for
% S5RPRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRRP8_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP8_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP8_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP8_gravloadJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t38 = sin(qJ(3));
t39 = sin(qJ(1));
t40 = cos(qJ(3));
t41 = cos(qJ(1));
t13 = -t39 * t38 - t41 * t40;
t25 = sin(qJ(4));
t26 = cos(qJ(4));
t31 = t26 * pkin(4) + t25 * qJ(5);
t51 = t31 * t13;
t14 = t41 * t38 - t39 * t40;
t50 = t31 * t14;
t32 = g(1) * t14 - g(2) * t13;
t49 = t32 * t25;
t48 = t32 * t26;
t34 = -t14 * pkin(3) - t13 * pkin(7);
t35 = t13 * pkin(3) - t14 * pkin(7);
t8 = g(1) * t13 + g(2) * t14;
t37 = t41 * pkin(1) + t39 * qJ(2);
t36 = t41 * pkin(2) + t37;
t33 = -t39 * pkin(1) + t41 * qJ(2);
t29 = t36 - t35;
t28 = -t39 * pkin(2) + t33;
t27 = t28 - t34;
t16 = g(1) * t41 + g(2) * t39;
t15 = g(1) * t39 - g(2) * t41;
t2 = -g(3) * t25 - t8 * t26;
t1 = g(3) * t26 - t8 * t25;
t3 = [0, 0, 0, 0, 0, 0, t15, t16, 0, 0, 0, 0, 0, 0, 0, 0, t15, 0, -t16, -g(1) * t33 - g(2) * t37, 0, 0, 0, 0, 0, 0, -t32, t8, 0, -g(1) * t28 - g(2) * t36, 0, 0, 0, 0, 0, 0, -t48, t49, -t8, -g(1) * t27 - g(2) * t29, 0, 0, 0, 0, 0, 0, -t48, -t8, -t49, -g(1) * (t27 + t50) - g(2) * (t29 - t51); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, -t8, 0, 0, 0, 0, 0, 0, 0, 0, t48, -t49, t8, -g(1) * t34 - g(2) * t35, 0, 0, 0, 0, 0, 0, t48, t8, t49, -g(1) * (t34 - t50) - g(2) * (t35 + t51); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t2, g(3) * t31 - t8 * (pkin(4) * t25 - qJ(5) * t26); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t3;
