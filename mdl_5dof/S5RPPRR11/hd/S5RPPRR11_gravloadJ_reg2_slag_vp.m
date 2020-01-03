% Calculate inertial parameters regressor of gravitation load for
% S5RPPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPPRR11_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR11_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR11_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR11_gravloadJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t15 = sin(qJ(4));
t18 = cos(qJ(4));
t22 = -t15 * pkin(4) + t18 * pkin(7);
t16 = sin(qJ(1));
t19 = cos(qJ(1));
t8 = g(1) * t19 + g(2) * t16;
t38 = -g(3) * t15 + t8 * t18;
t34 = g(3) * t18;
t31 = t19 * pkin(6);
t30 = t19 * pkin(1) + t16 * qJ(2);
t14 = sin(qJ(5));
t29 = t16 * t14;
t17 = cos(qJ(5));
t28 = t16 * t17;
t27 = t19 * t14;
t26 = t19 * t17;
t25 = -pkin(1) - qJ(3);
t24 = t19 * qJ(3) + t30;
t7 = g(1) * t16 - g(2) * t19;
t11 = t19 * qJ(2);
t21 = t25 * t16 + t11;
t6 = t7 * t18;
t5 = t15 * t26 - t29;
t4 = -t15 * t27 - t28;
t3 = -t15 * t28 - t27;
t2 = t15 * t29 - t26;
t1 = t8 * t15 + t34;
t9 = [0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, -t8, -g(1) * (-t16 * pkin(1) + t11) - g(2) * t30, 0, 0, 0, 0, 0, 0, 0, -t8, t7, -g(1) * t21 - g(2) * t24, 0, 0, 0, 0, 0, 0, t7 * t15, t6, t8, -g(1) * (t21 - t31) - g(2) * (-t16 * pkin(6) + t24), 0, 0, 0, 0, 0, 0, -g(1) * t3 - g(2) * t5, -g(1) * t2 - g(2) * t4, -t6, -g(1) * (t11 - t31) - g(2) * (-t22 * t19 + t24) + (-g(1) * (t22 + t25) + g(2) * pkin(6)) * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, t1, 0, 0, 0, 0, 0, 0, 0, 0, -t38 * t17, t38 * t14, -t1, -g(3) * t22 - t8 * (pkin(4) * t18 + pkin(7) * t15); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t4 + g(2) * t2 + t14 * t34, g(1) * t5 - g(2) * t3 + t17 * t34, 0, 0;];
taug_reg = t9;
