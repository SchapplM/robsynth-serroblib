% Calculate inertial parameters regressor of gravitation load for
% S5RRPPR4
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
% Datum: 2019-12-31 19:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPPR4_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR4_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR4_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR4_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t23 = sin(qJ(1));
t39 = t23 * pkin(1);
t35 = qJ(1) + qJ(2);
t20 = sin(t35);
t31 = cos(t35);
t38 = t31 * pkin(2) + t20 * qJ(3);
t37 = cos(pkin(8));
t36 = sin(pkin(8));
t17 = t31 * pkin(3);
t34 = t17 + t38;
t25 = cos(qJ(1));
t21 = t25 * pkin(1);
t33 = t21 + t38;
t16 = t31 * qJ(3);
t32 = -t20 * pkin(2) + t16;
t8 = -t20 * t36 - t31 * t37;
t9 = -t20 * t37 + t31 * t36;
t30 = g(1) * t9 - g(2) * t8;
t5 = g(1) * t8 + g(2) * t9;
t29 = t16 + (-pkin(2) - pkin(3)) * t20;
t28 = g(1) * t23 - g(2) * t25;
t27 = -t8 * pkin(4) + t9 * pkin(7) + t34;
t26 = t9 * pkin(4) + t8 * pkin(7) + t29;
t24 = cos(qJ(5));
t22 = sin(qJ(5));
t11 = g(1) * t31 + g(2) * t20;
t10 = g(1) * t20 - g(2) * t31;
t2 = t30 * t24;
t1 = t30 * t22;
t3 = [0, 0, 0, 0, 0, 0, t28, g(1) * t25 + g(2) * t23, 0, 0, 0, 0, 0, 0, 0, 0, t10, t11, 0, t28 * pkin(1), 0, 0, 0, 0, 0, 0, t10, 0, -t11, -g(1) * (t32 - t39) - g(2) * t33, 0, 0, 0, 0, 0, 0, -t30, t5, 0, -g(1) * (t29 - t39) - g(2) * (t17 + t33), 0, 0, 0, 0, 0, 0, -t2, t1, -t5, -g(1) * (t26 - t39) - g(2) * (t21 + t27); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, t11, 0, 0, 0, 0, 0, 0, 0, 0, t10, 0, -t11, -g(1) * t32 - g(2) * t38, 0, 0, 0, 0, 0, 0, -t30, t5, 0, -g(1) * t29 - g(2) * t34, 0, 0, 0, 0, 0, 0, -t2, t1, -t5, -g(1) * t26 - g(2) * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3) * t24 - t5 * t22, -g(3) * t22 - t5 * t24, 0, 0;];
taug_reg = t3;
