% Calculate minimal parameter regressor of gravitation load for
% S4PRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% taug_reg [4x18]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4PRRR6_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR6_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR6_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR6_gravloadJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t11 = cos(qJ(2));
t6 = sin(pkin(7));
t7 = cos(pkin(7));
t14 = g(1) * t7 + g(2) * t6;
t9 = sin(qJ(2));
t13 = -g(3) * t11 + t14 * t9;
t20 = g(3) * t9;
t18 = t11 * t6;
t17 = t11 * t7;
t8 = sin(qJ(3));
t16 = t11 * t8;
t10 = cos(qJ(3));
t15 = t10 * t11;
t5 = qJ(3) + qJ(4);
t4 = cos(t5);
t3 = sin(t5);
t2 = -g(1) * (-t4 * t17 - t6 * t3) - g(2) * (-t4 * t18 + t7 * t3) + t4 * t20;
t1 = -g(1) * (-t3 * t17 + t6 * t4) - g(2) * (-t3 * t18 - t7 * t4) + t3 * t20;
t12 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, t13, t14 * t11 + t20, 0, 0, 0, 0, 0, t13 * t10, -t13 * t8, 0, 0, 0, 0, 0, t13 * t4, -t13 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t6 * t10 - t7 * t16) - g(2) * (-t7 * t10 - t6 * t16) + t8 * t20, -g(1) * (-t7 * t15 - t6 * t8) - g(2) * (-t6 * t15 + t7 * t8) + t10 * t20, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg = t12;
