% Calculate inertial parameters regressor of gravitation load for
% S4RRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,theta3]';
% 
% Output:
% taug_reg [4x(4*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RRPP3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP3_gravloadJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP3_gravloadJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t16 = sin(qJ(1));
t18 = cos(qJ(1));
t6 = g(1) * t18 + g(2) * t16;
t5 = g(1) * t16 - g(2) * t18;
t13 = qJ(2) + pkin(6);
t10 = cos(t13);
t9 = sin(t13);
t21 = t10 * pkin(3) + t9 * qJ(4);
t15 = sin(qJ(2));
t17 = cos(qJ(2));
t19 = -g(3) * t17 + t6 * t15;
t14 = -qJ(3) - pkin(5);
t11 = t17 * pkin(2);
t8 = t11 + pkin(1);
t7 = t18 * t8;
t4 = t5 * t10;
t3 = t5 * t9;
t2 = g(3) * t9 + t6 * t10;
t1 = -g(3) * t10 + t6 * t9;
t12 = [0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, t5 * t17, -t5 * t15, -t6, -g(1) * (-t16 * pkin(1) + t18 * pkin(5)) - g(2) * (t18 * pkin(1) + t16 * pkin(5)), 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * (-t18 * t14 - t16 * t8) - g(2) * (-t16 * t14 + t7), 0, 0, 0, 0, 0, 0, t4, -t6, t3, -g(2) * t7 + (g(1) * t14 - g(2) * t21) * t18 + (-g(1) * (-t21 - t8) + g(2) * t14) * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, g(3) * t15 + t6 * t17, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t19 * pkin(2), 0, 0, 0, 0, 0, 0, t1, 0, -t2, -g(3) * (t11 + t21) + t6 * (pkin(2) * t15 + pkin(3) * t9 - qJ(4) * t10); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t12;
