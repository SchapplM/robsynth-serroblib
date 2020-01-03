% Calculate inertial parameters regressor of gravitation load for
% S4RPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% 
% Output:
% taug_reg [4x(4*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RPRP5_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP5_gravloadJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP5_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP5_gravloadJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t16 = sin(qJ(1));
t17 = cos(qJ(1));
t7 = g(1) * t17 + g(2) * t16;
t6 = g(1) * t16 - g(2) * t17;
t12 = pkin(6) + qJ(3);
t10 = cos(t12);
t9 = sin(t12);
t18 = t10 * pkin(3) + t9 * qJ(4);
t15 = -pkin(5) - qJ(2);
t14 = cos(pkin(6));
t8 = t14 * pkin(2) + pkin(1);
t5 = t17 * t8;
t4 = t6 * t10;
t3 = t6 * t9;
t2 = g(3) * t9 + t7 * t10;
t1 = -g(3) * t10 + t7 * t9;
t11 = [0, 0, 0, 0, 0, 0, t6, t7, 0, 0, 0, 0, 0, 0, 0, 0, t6 * t14, -t6 * sin(pkin(6)), -t7, -g(1) * (-t16 * pkin(1) + t17 * qJ(2)) - g(2) * (t17 * pkin(1) + t16 * qJ(2)), 0, 0, 0, 0, 0, 0, t4, -t3, -t7, -g(1) * (-t17 * t15 - t16 * t8) - g(2) * (-t16 * t15 + t5), 0, 0, 0, 0, 0, 0, t4, -t7, t3, -g(2) * t5 + (g(1) * t15 - g(2) * t18) * t17 + (-g(1) * (-t18 - t8) + g(2) * t15) * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t2, -g(3) * t18 + t7 * (pkin(3) * t9 - qJ(4) * t10); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t11;
