% Calculate inertial parameters regressor of gravitation load for
% S4PRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% 
% Output:
% taug_reg [4x(4*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4PRRP5_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP5_gravloadJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP5_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP5_gravloadJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t8 = sin(pkin(6));
t9 = cos(pkin(6));
t17 = g(1) * t9 + g(2) * t8;
t12 = sin(qJ(2));
t14 = cos(qJ(2));
t5 = -g(3) * t14 + t17 * t12;
t21 = g(3) * t12;
t11 = sin(qJ(3));
t19 = t11 * t14;
t13 = cos(qJ(3));
t18 = t13 * t14;
t1 = -g(1) * (t8 * t13 - t9 * t19) - g(2) * (-t9 * t13 - t8 * t19) + t11 * t21;
t10 = -qJ(4) - pkin(5);
t7 = t13 * pkin(3) + pkin(2);
t6 = t17 * t14 + t21;
t4 = t5 * t13;
t3 = t5 * t11;
t2 = -g(1) * (-t8 * t11 - t9 * t18) - g(2) * (t9 * t11 - t8 * t18) + t13 * t21;
t15 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(3) * (t14 * pkin(2) + t12 * pkin(5)) + t17 * (pkin(2) * t12 - pkin(5) * t14), 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(3) * (-t12 * t10 + t14 * t7) + t17 * (t10 * t14 + t12 * t7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5;];
taug_reg = t15;
