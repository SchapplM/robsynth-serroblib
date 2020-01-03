% Calculate inertial parameters regressor of gravitation load for
% S4RRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% 
% Output:
% taug_reg [4x(4*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RRPR5_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR5_gravloadJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR5_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR5_gravloadJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t13 = qJ(1) + qJ(2);
t10 = sin(t13);
t11 = cos(t13);
t23 = t11 * pkin(2) + t10 * qJ(3);
t15 = sin(qJ(1));
t22 = t15 * pkin(1);
t17 = cos(qJ(1));
t21 = t17 * pkin(1) + t23;
t6 = t11 * qJ(3);
t20 = -t10 * pkin(2) + t6;
t19 = t6 + (-pkin(2) - pkin(6)) * t10;
t4 = g(1) * t11 + g(2) * t10;
t3 = g(1) * t10 - g(2) * t11;
t18 = g(1) * t15 - g(2) * t17;
t16 = cos(qJ(4));
t14 = sin(qJ(4));
t7 = t11 * pkin(6);
t2 = t4 * t16;
t1 = t4 * t14;
t5 = [0, 0, 0, 0, 0, 0, t18, g(1) * t17 + g(2) * t15, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, t18 * pkin(1), 0, 0, 0, 0, 0, 0, 0, -t3, -t4, -g(1) * (t20 - t22) - g(2) * t21, 0, 0, 0, 0, 0, 0, -t1, -t2, t3, -g(1) * (t19 - t22) - g(2) * (t7 + t21); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, -t4, -g(1) * t20 - g(2) * t23, 0, 0, 0, 0, 0, 0, -t1, -t2, t3, -g(1) * t19 - g(2) * (t7 + t23); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3) * t14 - t3 * t16, g(3) * t16 + t3 * t14, 0, 0;];
taug_reg = t5;
