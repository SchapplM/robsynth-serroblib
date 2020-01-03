% Calculate inertial parameters regressor of gravitation load for
% S5RPPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta3,theta4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPPPR5_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR5_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR5_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR5_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t31 = cos(qJ(1));
t30 = sin(qJ(1));
t29 = t31 * pkin(1) + t30 * qJ(2);
t28 = cos(pkin(7));
t27 = sin(pkin(7));
t26 = t31 * pkin(2) + t29;
t3 = -t30 * t27 - t31 * t28;
t4 = t31 * t27 - t30 * t28;
t25 = g(1) * t4 - g(2) * t3;
t24 = g(1) * t3 + g(2) * t4;
t23 = -t30 * pkin(1) + t31 * qJ(2);
t22 = -t30 * pkin(2) + t23;
t21 = -pkin(6) - qJ(4);
t20 = cos(pkin(8));
t18 = pkin(8) + qJ(5);
t12 = cos(t18);
t11 = sin(t18);
t10 = t20 * pkin(4) + pkin(3);
t6 = g(1) * t31 + g(2) * t30;
t5 = g(1) * t30 - g(2) * t31;
t1 = [0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, t5, 0, -t6, -g(1) * t23 - g(2) * t29, 0, 0, 0, 0, 0, 0, -t25, t24, 0, -g(1) * t22 - g(2) * t26, 0, 0, 0, 0, 0, 0, -t25 * t20, t25 * sin(pkin(8)), -t24, -g(1) * (t4 * pkin(3) + t3 * qJ(4) + t22) - g(2) * (-t3 * pkin(3) + t4 * qJ(4) + t26), 0, 0, 0, 0, 0, 0, -t25 * t12, t25 * t11, -t24, -g(1) * (t4 * t10 - t3 * t21 + t22) - g(2) * (-t3 * t10 - t4 * t21 + t26); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3) * t12 - t24 * t11, -g(3) * t11 - t24 * t12, 0, 0;];
taug_reg = t1;
