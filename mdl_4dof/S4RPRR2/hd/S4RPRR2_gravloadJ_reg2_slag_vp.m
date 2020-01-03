% Calculate inertial parameters regressor of gravitation load for
% S4RPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% taug_reg [4x(4*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RPRR2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR2_gravloadJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR2_gravloadJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t15 = qJ(1) + pkin(7);
t13 = qJ(3) + t15;
t10 = cos(t13);
t9 = sin(t13);
t24 = t10 * pkin(3) + t9 * pkin(6);
t12 = cos(t15);
t19 = cos(qJ(1));
t23 = t19 * pkin(1) + pkin(2) * t12;
t22 = -t9 * pkin(3) + t10 * pkin(6);
t3 = g(1) * t9 - g(2) * t10;
t4 = g(1) * t10 + g(2) * t9;
t11 = sin(t15);
t17 = sin(qJ(1));
t21 = -t17 * pkin(1) - pkin(2) * t11;
t20 = g(1) * t17 - g(2) * t19;
t18 = cos(qJ(4));
t16 = sin(qJ(4));
t2 = t3 * t18;
t1 = t3 * t16;
t5 = [0, 0, 0, 0, 0, 0, t20, g(1) * t19 + g(2) * t17, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t11 - g(2) * t12, g(1) * t12 + g(2) * t11, 0, t20 * pkin(1), 0, 0, 0, 0, 0, 0, t3, t4, 0, -g(1) * t21 - g(2) * t23, 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * (t21 + t22) - g(2) * (t23 + t24); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * t22 - g(2) * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t18 + t4 * t16, g(3) * t16 + t4 * t18, 0, 0;];
taug_reg = t5;
