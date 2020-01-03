% Calculate inertial parameters regressor of gravitation load for
% S4PRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
% 
% Output:
% taug_reg [4x(4*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4PRPR6_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR6_gravloadJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR6_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR6_gravloadJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t11 = cos(pkin(6));
t9 = sin(pkin(6));
t18 = g(1) * t11 + g(2) * t9;
t13 = sin(qJ(2));
t14 = cos(qJ(2));
t1 = -g(3) * t14 + t18 * t13;
t22 = g(3) * t13;
t20 = t14 * t9;
t19 = t11 * t14;
t12 = -pkin(5) - qJ(3);
t10 = cos(pkin(7));
t7 = pkin(7) + qJ(4);
t5 = cos(t7);
t4 = sin(t7);
t3 = t10 * pkin(3) + pkin(2);
t2 = t18 * t14 + t22;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1 * t10, -t1 * sin(pkin(7)), -t2, -g(3) * (t14 * pkin(2) + t13 * qJ(3)) + t18 * (pkin(2) * t13 - qJ(3) * t14), 0, 0, 0, 0, 0, 0, t1 * t5, -t1 * t4, -t2, -g(3) * (-t13 * t12 + t14 * t3) + t18 * (t12 * t14 + t13 * t3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t4 * t19 + t9 * t5) - g(2) * (-t11 * t5 - t4 * t20) + t4 * t22, -g(1) * (-t5 * t19 - t9 * t4) - g(2) * (t11 * t4 - t5 * t20) + t5 * t22, 0, 0;];
taug_reg = t6;
