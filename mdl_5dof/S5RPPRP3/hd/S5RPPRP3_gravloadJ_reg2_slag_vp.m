% Calculate inertial parameters regressor of gravitation load for
% S5RPPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPPRP3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP3_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP3_gravloadJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t16 = sin(qJ(4));
t23 = pkin(4) * t16;
t14 = qJ(1) + pkin(7);
t11 = sin(t14);
t12 = cos(t14);
t19 = cos(qJ(1));
t22 = t19 * pkin(1) + t12 * pkin(2) + t11 * qJ(3);
t17 = sin(qJ(1));
t21 = -t17 * pkin(1) + t12 * qJ(3);
t6 = g(1) * t12 + g(2) * t11;
t5 = g(1) * t11 - g(2) * t12;
t20 = g(1) * t17 - g(2) * t19;
t18 = cos(qJ(4));
t2 = g(3) * t16 - t5 * t18;
t15 = -qJ(5) - pkin(6);
t4 = t6 * t18;
t3 = t6 * t16;
t1 = g(3) * t18 + t5 * t16;
t7 = [0, 0, 0, 0, 0, 0, t20, g(1) * t19 + g(2) * t17, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, t20 * pkin(1), 0, 0, 0, 0, 0, 0, 0, -t5, -t6, -g(1) * (-t11 * pkin(2) + t21) - g(2) * t22, 0, 0, 0, 0, 0, 0, -t3, -t4, t5, -g(1) * ((-pkin(2) - pkin(6)) * t11 + t21) - g(2) * (t12 * pkin(6) + t22), 0, 0, 0, 0, 0, 0, -t3, -t4, t5, -g(1) * (t12 * t23 + (-pkin(2) + t15) * t11 + t21) - g(2) * (t11 * t23 - t12 * t15 + t22); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6;];
taug_reg = t7;
