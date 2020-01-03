% Calculate inertial parameters regressor of gravitation load for
% S5RPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRRR4_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR4_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t25 = qJ(1) + pkin(9);
t22 = qJ(3) + t25;
t19 = qJ(4) + t22;
t13 = sin(t19);
t14 = cos(t19);
t37 = t14 * pkin(4) + t13 * pkin(8);
t20 = sin(t25);
t27 = sin(qJ(1));
t36 = t27 * pkin(1) + pkin(2) * t20;
t21 = cos(t25);
t29 = cos(qJ(1));
t35 = t29 * pkin(1) + pkin(2) * t21;
t18 = cos(t22);
t12 = pkin(3) * t18;
t34 = t12 + t35;
t33 = t13 * pkin(4) - t14 * pkin(8);
t17 = sin(t22);
t11 = pkin(3) * t17;
t32 = t11 + t33;
t31 = g(2) * t14 + g(3) * t13;
t3 = g(2) * t13 - g(3) * t14;
t6 = -g(2) * t18 - g(3) * t17;
t30 = -g(2) * t29 - g(3) * t27;
t28 = cos(qJ(5));
t26 = sin(qJ(5));
t5 = g(2) * t17 - g(3) * t18;
t2 = t31 * t28;
t1 = t31 * t26;
t4 = [0, 0, 0, 0, 0, 0, t30, g(2) * t27 - g(3) * t29, 0, 0, 0, 0, 0, 0, 0, 0, -g(2) * t21 - g(3) * t20, g(2) * t20 - g(3) * t21, 0, t30 * pkin(1), 0, 0, 0, 0, 0, 0, t6, t5, 0, -g(2) * t35 - g(3) * t36, 0, 0, 0, 0, 0, 0, -t31, t3, 0, -g(2) * t34 - g(3) * (t11 + t36), 0, 0, 0, 0, 0, 0, -t2, t1, -t3, -g(2) * (t34 + t37) - g(3) * (t32 + t36); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, t5, 0, 0, 0, 0, 0, 0, 0, 0, -t31, t3, 0, t6 * pkin(3), 0, 0, 0, 0, 0, 0, -t2, t1, -t3, -g(2) * (t12 + t37) - g(3) * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, t3, 0, 0, 0, 0, 0, 0, 0, 0, -t2, t1, -t3, -g(2) * t37 - g(3) * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t28 + t3 * t26, g(1) * t26 + t3 * t28, 0, 0;];
taug_reg = t4;
